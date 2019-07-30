<template>
  <v-content>
    <v-container fluid fill-height>
      <v-layout align-center justify-center>
        <v-flex xs12 sm8 md4>
          <v-card class="elevation-12">
            <v-toolbar dark color="primary">
              <v-toolbar-title>Create {{appName}} Account</v-toolbar-title>
              <v-spacer></v-spacer>
            </v-toolbar>
            <v-card-text>
              <v-form
                @keyup.enter="submit"
                v-model="valid"
                ref="form"
                @submit.prevent=""
                lazy-validation
              >
                <v-text-field
                  @keyup.enter="submit"
                  type="email"
                  v-model="email"
                  prepend-icon="mdi-account"
                  name="email"
                  label="Email"
                  data-vv-as="email"
                  v-validate="'required|email'"
                  :error-messages="errors.first('email')"
                ></v-text-field>
                <v-text-field
                  @keyup.enter="submit"
                  type="text"
                  label="Full Name"
                  prepend-icon="mdi-account"
                  v-model="fullName"
                ></v-text-field>
                <v-text-field
                  @keyup.enter="submit"
                  type="password"
                  ref="password"
                  label="Password"
                  data-vv-name="password"
                  data-vv-delay="100"
                  data-vv-rules="required"
                  v-validate="'required'"
                  v-model="password1"
                  prepend-icon="mdi-lock"
                  :error-messages="errors.first('password')"
                  required
                ></v-text-field>
                <v-text-field
                  @keyup.enter="submit"
                  type="password"
                  label="Confirm Password"
                  data-vv-name="password_confirmation"
                  data-vv-delay="100"
                  data-vv-rules="required|confirmed:$password"
                  data-vv-as="password"
                  v-validate="'required|confirmed:password'"
                  v-model="password2"
                  prepend-icon="mdi-lock"
                  :error-messages="errors.first('password_confirmation')"
                  required
                ></v-text-field>
              </v-form>
              <v-flex class="caption text-right">
                <router-link to="/login">Already have an account?</router-link>
              </v-flex>
            </v-card-text>
            <v-card-actions>
              <v-spacer></v-spacer>
              <v-btn @click.prevent="submit">Sign Up</v-btn>
            </v-card-actions>
          </v-card>
        </v-flex>
      </v-layout>
    </v-container>
  </v-content>
</template>

<script lang="ts">
  import { appName } from '@/env';
  import { mainModule } from '@/modules/main';
  import { userModule } from '@/modules/user';
  import { Component, Vue } from 'vue-property-decorator';

  @Component
  export default class SignUp extends Vue {
    readonly mainContext = mainModule.context(this.$store);
    readonly userContext = userModule.context(this.$store);

    valid = true;
    email = '';
    fullName = '';
    password1 = '';
    password2 = '';
    appName = appName;

    reset() {
      this.email = '';
      this.fullName = '';
      this.password1 = '';
      this.password2 = '';
      this.$validator.reset();
    }

    async submit() {
      if (await this.$validator.validateAll()) {
        const userExist = await this.checkUserExists();
        if (!userExist) {
          await this.userContext.actions.signUp({ email: this.email, password: this.password1 });
        }
      }
    }

    private async checkUserExists() {
      const userExist = await this.userContext.actions.checkUserExists(this.email);
      if (userExist) {
        this.mainContext.mutations.addNotification({
          content: 'User with this email already exists',
          color: 'warning',
        });
      }
      return userExist;
    }
  }
</script>

<style>
</style>
