const axiosConfig = {
    timeout: 30000,
    withCredentials: false
};


import Vue from 'vue'
import App from './App.vue'
import axios from 'axios'
import VueFlashMessage from 'vue-flash-message';
require('vue-flash-message/dist/vue-flash-message.min.css')

Vue.config.productionTip = false
Vue.prototype.$axios = axios.create(axiosConfig);
Vue.use(VueFlashMessage);

new Vue({
  render: function (h) { return h(App) },
}).$mount('#app')
