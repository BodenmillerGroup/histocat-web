<template>
  <v-container fluid>
    <v-card class="ma-4 pa-4">
      <v-card-title primary-title>
        <div class="headline primary--text">Edit User</div>
      </v-card-title>
      <v-card-text>
        <template>
          <div class="my-4">
            <div class="subtitle-1 secondary--text text--lighten-2">Username</div>
            <div
              class="title primary--text text--darken-2"
              v-if="user"
            >{{user.email}}
            </div>
            <div
              class="title primary--text text--darken-2"
              v-else
            >-----
            </div>
          </div>
          <v-form
            v-model="valid"
            ref="form"
            lazy-validation
          >
            <v-text-field
              label="Full Name"
              v-model="fullName"
              required
            ></v-text-field>
            <v-text-field
              label="E-mail"
              type="email"
              v-model="email"
              v-validate="'required|email'"
              data-vv-name="email"
              :error-messages="errors.collect('email')"
              required
            ></v-text-field>
            <div class="subtitle-1 secondary--text text--lighten-2">User is superuser <span v-if="isSuperuser">(currently is a superuser)</span><span
              v-else>(currently is not a superuser)</span></div>
            <v-checkbox
              label="Is Superuser"
              v-model="isSuperuser"
            ></v-checkbox>
            <div class="subtitle-1 secondary--text text--lighten-2">User is active <span v-if="isActive">(currently active)</span><span
              v-else>(currently not active)</span></div>
            <v-checkbox
              label="Is Active"
              v-model="isActive"
            ></v-checkbox>
            <v-layout align-center>
              <v-flex shrink>
                <v-checkbox
                  v-model="setPassword"
                  class="mr-2"
                ></v-checkbox>
              </v-flex>
              <v-flex>
                <v-text-field
                  :disabled="!setPassword"
                  type="password"
                  ref="password"
                  label="Set Password"
                  data-vv-name="password"
                  data-vv-delay="100"
                  v-validate="{required: setPassword}"
                  v-model="password1"
                  :error-messages="errors.first('password')"
                >
                </v-text-field>
                <v-text-field
                  v-show="setPassword"
                  type="password"
                  label="Confirm Password"
                  data-vv-name="password_confirmation"
                  data-vv-delay="100"
                  data-vv-as="password"
                  v-validate="{required: setPassword, confirmed: 'password'}"
                  v-model="password2"
                  :error-messages="errors.first('password_confirmation')"
                >
                </v-text-field>
              </v-flex>
            </v-layout>
          </v-form>
        </template>
      </v-card-text>
      <v-card-actions>
        <v-spacer></v-spacer>
        <v-btn @click="cancel">Cancel</v-btn>
        <v-btn @click="reset">Reset</v-btn>
        <v-btn
          @click="submit"
          :disabled="!valid"
        >
          Save
        </v-btn>
      </v-card-actions>
    </v-card>
  </v-container>
</template>

<script lang="ts">
  import { userModule } from '@/modules/user';
  import { IUserProfileUpdate } from '@/modules/user/models';
  import { Component, Vue } from 'vue-property-decorator';

  @Component
  export default class EditUser extends Vue {
    readonly userContext = userModule.context(this.$store);

    valid = true;
    fullName: string = '';
    email: string = '';
    isActive: boolean = true;
    isSuperuser: boolean = false;
    setPassword = false;
    password1: string = '';
    password2: string = '';

    async mounted() {
      await this.userContext.actions.getUsers();
      this.reset();
    }

    reset() {
      this.setPassword = false;
      this.password1 = '';
      this.password2 = '';
      this.$validator.reset();
      if (this.user) {
        this.fullName = this.user.full_name;
        this.email = this.user.email;
        this.isActive = this.user.is_active;
        this.isSuperuser = this.user.is_superuser;
      }
    }

    cancel() {
      this.$router.back();
    }

    async submit() {
      if (await this.$validator.validateAll()) {
        const data: IUserProfileUpdate = {};
        if (this.fullName) {
          data.full_name = this.fullName;
        }
        if (this.email) {
          data.email = this.email;
        }
        data.is_active = this.isActive;
        data.is_superuser = this.isSuperuser;
        if (this.setPassword) {
          data.password = this.password1;
        }
        await this.userContext.actions.updateUser({ id: this.user!.id, user: data });
        this.$router.push('/main/admin/users');
      }
    }

    get user() {
      return this.userContext.getters.adminOneUser(+this.$router.currentRoute.params.id);
    }
  }
</script>
